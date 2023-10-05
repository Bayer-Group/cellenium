import { atom, selector } from 'recoil';
import * as aq from 'arquero';
import { Struct } from 'arquero/dist/types/table/transformable';
import { Omics, SamplesAnnotationTable, SamplesProjectionTable, Study, StudyOmicsTable } from './model';
import { apolloClient } from './client';
import {
  AllGenesDocument,
  AllGenesQuery,
  AllGenesQueryVariables,
  CellOAnnotationGroupIdDocument,
  CellOAnnotationGroupIdQuery,
  CellOAnnotationGroupIdQueryVariables,
  OmicsGeneFragment,
  StudyAnnotationFrontendValue,
  StudyBasicsDocument,
  StudyBasicsFragment,
  StudyBasicsQuery,
  StudyBasicsQueryVariables,
} from './generated/types';

export const studyIdState = atom<number>({
  key: 'studyId',
  default: undefined,
});

function buildSampleAnnotationTable(s: StudyBasicsFragment) {
  const samplesTable = aq
    .from(s.studySampleAnnotationSubsamplingList)
    .select(aq.not(['__typename']))
    .unroll('studySampleIds')
    .select({
      studySampleIds: 'studySampleId',
      annotationValueId: 'annotationValueId',
    });
  // samplesTable.print();
  const annotationGroupsValuesTable = aq
    .from(s.annotationGroupsList)
    .unroll('annotationValuesList')
    .derive({
      // @ts-ignore
      annotationValueId: (d) => d.annotationValuesList.annotationValueId,
    })
    .select('annotationGroupId', 'annotationValueId');

  // annotationGroupsValuesTable.print();
  const annotatedSamplesTable = samplesTable.join(annotationGroupsValuesTable, 'annotationValueId').reify();
  return SamplesAnnotationTable.definedTable(annotatedSamplesTable);
}

function buildOmicsTable(d: { displaySymbol: string[]; displayName: string[]; omicsType: string[]; omicsId: number[] }) {
  const retTable = aq.table({
    value: d.displaySymbol,
    displaySymbol: d.displaySymbol,
    displayName: d.displayName,
    omicsType: d.omicsType,
    omicsId: d.omicsId,
  });
  return StudyOmicsTable.definedTable(retTable);
}

function buildSampleProjectionTable(d: { studySampleId: number[]; projection: number[] }) {
  return SamplesProjectionTable.definedTable(
    aq.table({
      studySampleId: d.studySampleId,
      projectionX: Array.from(Array(d.projection ? d.projection.length / 2 : 0).keys()).map((i) => d.projection[i * 2]),
      projectionY: Array.from(Array(d.projection ? d.projection.length / 2 : 0).keys()).map((i) => d.projection[i * 2 + 1]),
    }),
  );
}

export const studyState = selector<Study | undefined>({
  key: 'studyState',
  get: async ({ get }) => {
    const studyId = get(studyIdState);
    const studyReloadHelper = get(studyReloadHelperState);
    if (studyId && studyReloadHelper) {
      const responsePromise = apolloClient.query<StudyBasicsQuery, StudyBasicsQueryVariables>({
        query: StudyBasicsDocument,
        variables: {
          studyId,
        },
        fetchPolicy: 'network-only',
      });
      // could do multiple queries in parallel ... but maybe not needed
      const response = await responsePromise;

      if (response?.data?.study) {
        if (response.data.study.studySampleProjectionSubsamplingTransposedList[0].projection.length === 0) {
          throw Error('no projection data');
        }
        if (response.data.study.studySampleAnnotationSubsamplingList.length === 0) {
          throw Error('no study annotations');
        }
        if (response.data.study.studyOmicsTransposedList.length === 0) {
          throw Error('no genes');
        }
        const studyOmicsTable = buildOmicsTable(response.data.study.studyOmicsTransposedList[0]);
        const omicsTypes = studyOmicsTable.rollup({ omicsTypes: aq.op.array_agg_distinct('omicsType') }).array('omicsTypes')[0];
        const s: Study = {
          ...response.data.study,
          samplesProjectionTables: new Map(
            response.data.study.studySampleProjectionSubsamplingTransposedList.map((tl) => [
              tl.modality ? `${tl.modality}:${tl.projectionType}` : tl.projectionType,
              buildSampleProjectionTable(tl),
            ]),
          ),
          samplesAnnotationTable: buildSampleAnnotationTable(response.data.study),
          studyOmicsTable,
          studyOmicsMap: new Map(studyOmicsTable.objects().map((o) => [(o as Omics).omicsId, o as Omics])),
          omicsTypes,
          annotationGroupMap: new Map(response.data.study.annotationGroupsList.map((g) => [g.annotationGroupId, g])),
          annotationValueMap: new Map<number, StudyAnnotationFrontendValue>(
            response.data.study.annotationGroupsList
              .map((g) => g.annotationValuesList)
              .flat(2)
              .map((v: Struct) => [v.annotationValueId, v as StudyAnnotationFrontendValue]),
          ),
        };
        // console.log('studyState set', s);
        return s;
      }
      return undefined;
    }
    return null;
  },
  // by default, recoil protects returned objects with immutability - but arquero's query builder pattern needs to write to the table state
  dangerouslyAllowMutability: true,
});

export const pageState = atom<string>({
  key: 'page',
  default: 'CellMarkerAnalysis',
});

export const allGenesState = selector<Map<number, OmicsGeneFragment> | undefined>({
  key: 'allGenesState',
  get: async () => {
    const allGenes = await apolloClient.query<AllGenesQuery, AllGenesQueryVariables>({
      query: AllGenesDocument,
      fetchPolicy: 'no-cache',
    });
    if (allGenes?.data) {
      return new Map(allGenes.data.omicsBasesList.map((o: OmicsGeneFragment) => [o.omicsId, o]));
    }
    throw Error('Error fetching genes');
  },
});

export const cellOAnnotationGroupIdState = selector<number | undefined>({
  key: 'cellOAnnotationGroupIdState',
  get: async () => {
    const annotationGroupIdData = await apolloClient.query<CellOAnnotationGroupIdQuery, CellOAnnotationGroupIdQueryVariables>({
      query: CellOAnnotationGroupIdDocument,
      fetchPolicy: 'no-cache',
    });
    if (annotationGroupIdData?.data) {
      return annotationGroupIdData.data.annotationGroupsList[0].annotationGroupId;
    }
    return 0;
  },
});

export const correlationOmicsIdState = atom<number>({
  key: 'coexpressionOmicsId',
  default: undefined,
});

export const userGeneStoreOpenState = atom<boolean>({
  key: 'usergenestoreopen',
  default: false,
});
export const userGeneStoreCounterColor = atom<string>({
  key: 'usergenestorecountercolor',
  default: 'blue',
});

export const selectedGenesState = atom<Omics[]>({
  key: 'selectedGenes',
  default: [],
});

export const userGenesState = atom<Omics[]>({
  key: 'userGenes',
  default: [],
});

export const celltypeDiscoveryGenesState = atom<(Omics | null)[]>({
  key: 'celltypeDiscoveryGenesState',
  default: [null, null],
});

export const celltypeDiscoveryCoexpressionSamplesState = atom<(number[] | null)[]>({
  key: 'celltypeDiscoveryCoexpressionSamplesState',
  default: [null],
});

export const selectedAnnotationFilterState = atom<number[]>({
  key: 'selectedAnnotationFilter',
  default: [],
});

export const selectedAnnotationState = atom<number>({
  key: 'selectedAnnotation',
  default: 0,
});

export const highlightAnnotationState = atom<number>({
  key: 'highlightAnnotation',
  default: 0,
});

export const annotationGroupIdState = atom<number | undefined>({
  key: 'annotationGroupId',
  default: undefined,
});

export const annotationSecondaryGroupIdState = atom<number | undefined>({
  key: 'annotationSecondaryGroupId',
  default: undefined,
});

export const studyReloadHelperState = atom<number>({
  key: 'studyReloadHelper',
  default: 1,
});

// export const studyLayerIdDefinedState = atom<number>({
//     key: "studyLayerIdDefined",
//     default: undefined
// });

export const studyLayerIdState = selector<number>({
  key: 'studyLayerId',
  get: ({ get }) => {
    const study = get(studyState);
    if (study) {
      // check if the selected layer (additional atom studyLayerIdDefinedState, unused) is valid or return default:
      return study.studyLayersList[0].studyLayerId;
    }
    return -1;
  },
});

export const selectedProjectionState = atom<string>({
  key: 'selectedProjection',
  default: '',
});
