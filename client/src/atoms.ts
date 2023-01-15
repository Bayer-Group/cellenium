import {atom, selector} from "recoil";
import {Omics, SamplesAnnotationTable, SamplesProjectionTable, Study, StudyOmicsTable} from "./model";
import {apolloClient} from "./index";
import {StudyBasicsDocument, StudyBasicsFragment, StudyBasicsQuery, StudyBasicsQueryVariables} from "./generated/types";
import * as aq from 'arquero';

export const selectedGenesState = atom<Omics[]>({
    key: "selectedGenes",
    default: []
})

export const userGenesState = atom<Omics[]>({
    key: "userGenes",
    default: []
})

export const selectedAnnotationFilterState = atom<number[]>({
    key: "selectedAnnotationFilter",
    default: []
})

export const selectedAnnotationState = atom<number | undefined>({
    key: "selectedAnnotation",
    default: undefined
})

export const highlightAnnotationState = atom<number>({
    key: "highlightAnnotation",
    default: 0
})

export const annotationGroupIdState = atom<number | undefined>({
    key: "annotationGroupId",
    default: undefined
})

export const studyIdState = atom<number>({
    key: "studyId",
    default: undefined
});

export const studyLayerIdDefinedState = atom<number>({
    key: "studyLayerIdDefined",
    default: undefined
});

export const studyLayerIdState = selector<number>({
    key: "studyLayerId",
    get: ({get}) => {
        const study = get(studyState);
        if (study) {
            // check if the selected layer (additional atom) is valid or return default:
            return study.studyLayersList[0].studyLayerId;
        }
        return -1;
    }
});

function buildSampleProjectionTable(d: { studySampleId: number[], projection: number[] }) {
    return SamplesProjectionTable.definedTable(aq.table({
        studySampleId: d.studySampleId,
        projectionX: Array.from(Array(d.projection?.length / 2).keys()).map(i => d.projection[i * 2]),
        projectionY: Array.from(Array(d.projection?.length / 2).keys()).map(i => d.projection[i * 2 + 1])
    }));
}

function buildSampleAnnotationTable(s: StudyBasicsFragment) {
    const samplesTable = aq.from(s.studySampleAnnotationSubsamplingList).select(aq.not(['__typename']))
        .unroll('studySampleIds')
        .select({studySampleIds: 'studySampleId', annotationValueId: 'annotationValueId'});
    // samplesTable.print();
    const annotationGroupsValuesTable = aq.from(s.annotationGroupsList)
        .unroll('annotationValuesList')
        // @ts-ignore
        .derive({annotationValueId: r => r.annotationValuesList.annotationValueId})
        .select('annotationGroupId', 'annotationValueId');
    // annotationGroupsValuesTable.print();
    const annotatedSamplesTable = samplesTable.join(annotationGroupsValuesTable, 'annotationValueId').reify();
    return SamplesAnnotationTable.definedTable(annotatedSamplesTable);
}

function buildOmicsTable(d: { displaySymbol: string[], displayName: string[], omicsType: string[], omicsId: number[] }) {
    const retTable = aq.table({
        value: d.displaySymbol,
        displaySymbol: d.displaySymbol,
        displayName: d.displayName,
        omicsType: d.omicsType,
        omicsId: d.omicsId
    })
    return StudyOmicsTable.definedTable(retTable)
}

export const studyState = selector<Study | undefined>({
    key: "studyState",
    get: async ({get}) => {
        const studyId = get(studyIdState);
        if (studyId) {
            const responsePromise = apolloClient.query<StudyBasicsQuery, StudyBasicsQueryVariables>({
                query: StudyBasicsDocument,
                variables: {
                    studyId,
                },
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

                // do some computations, e.g. generate arquero table of a received record list...
                const s: Study = {
                    ...response.data.study,
                    samplesProjectionTable: buildSampleProjectionTable(response.data.study.studySampleProjectionSubsamplingTransposedList[0]),
                    samplesAnnotationTable: buildSampleAnnotationTable(response.data.study),
                    studyOmicsTable: buildOmicsTable(response.data.study.studyOmicsTransposedList[0]),
                    annotationGroupMap: new Map(response.data.study.annotationGroupsList.map((g) => [g.annotationGroupId, g])),
                    annotationValueMap: new Map(response.data.study.annotationGroupsList.map((g) => g.annotationValuesList).flat(2).map((v: any) => [v.annotationValueId, v]))
                };
                return s;
            }
        }
    },
    // by default, recoil protects returned objects with immutability - but arquero's query builder pattern needs to write to the table state
    dangerouslyAllowMutability: true
});


export const pageState = atom<string>({
    key: "page",
    default: "CellMarkerAnalysis"
});
