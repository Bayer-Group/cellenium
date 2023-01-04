import {atom, selector} from "recoil";
import {Study} from "./model";
import {apolloClient} from "./index";
import {
    StudyBasicsDocument, StudyBasicsFragment,
    StudyBasicsQuery,
    StudyBasicsQueryVariables,
    StudySampleProjectionSubsamplingTransposed
} from "./generated/types";
import * as aq from 'arquero';
import ColumnTable from 'arquero/dist/types/table/column-table';

export const studyIdState = atom<number | undefined>({
    key: "studyId",
    default: undefined,
});

function buildSampleProjectionTable(d: { studySampleId: number[], projection: number[] }) {
    return aq.table({
        studySampleId: d.studySampleId,
        projectionX: Array.from(Array(d.projection?.length / 2).keys()).map(i => d.projection[i * 2]),
        projectionY: Array.from(Array(d.projection?.length / 2).keys()).map(i => d.projection[i * 2 + 1])
    });
}

function buildSampleAnnotationTable(s: StudyBasicsFragment) {
    const samplesTable = aq.from(s.studySampleAnnotationsList).select(aq.not(['__typename']))
        .unroll('studySampleIds')
        .select({studySampleIds: 'studySampleId', annotationValueId: 'annotationValueId'});
    // samplesTable.print();
    const annotationGroupsValuesTable = aq.from(s.studyAnnotationGroupUisList.map(e => e.annotationGroup))
        .unroll('annotationValuesList')
        // @ts-ignore
        .derive({annotationValueId: r => r.annotationValuesList.annotationValueId})
        .select('annotationGroupId', 'annotationValueId');
    // annotationGroupsValuesTable.print();
    const annotatedSamplesTable = samplesTable.join(annotationGroupsValuesTable, 'annotationValueId').reify();
    annotatedSamplesTable.print();
    return annotatedSamplesTable;
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
                // do some computations, e.g. generate arquero table of a received record list...
                const s: Study = {
                    ...response.data.study,
                    samplesProjectionTable: buildSampleProjectionTable(response.data.study.studySampleProjectionSubsamplingTransposedList[0]),
                    samplesAnnotationTable: buildSampleAnnotationTable(response.data.study)
                };
                return s;
            }
        }
    },
    // by default, recoil protects returned objects with immutability - but arquero's query builder pattern needs to write to the table state
    dangerouslyAllowMutability: true
});
