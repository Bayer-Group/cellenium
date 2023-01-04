import {atom, selector} from "recoil";
import {Study} from "./model";
import {apolloClient} from "./index";
import {
    StudyBasicsDocument,
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
        sampleId: d.studySampleId,
        projectionX: Array.from(Array(d.projection?.length / 2).keys()).map(i => d.projection[i * 2]),
        projectionY: Array.from(Array(d.projection?.length / 2).keys()).map(i => d.projection[i * 2 + 1])
    });
}

function dummyStudy(): Study {
    return {
        samplesProjectionTable: aq.table({})
    } as Study;
}

export const studyState = selector<Study>({
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
                    samplesProjectionTable: buildSampleProjectionTable(response.data.study.studySampleProjectionSubsamplingTransposedList[0])
                };
                return s;
            }
            return dummyStudy();
        } else {
            return dummyStudy();
        }
    }
});
