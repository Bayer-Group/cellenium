import {atom, selector} from "recoil";
import {Study} from "./model";
import {apolloClient} from "./index";
import {StudyBasicsDocument, StudyBasicsQuery, StudyBasicsQueryVariables} from "./generated/types";

export const studyIdState = atom<number | undefined>({
    key: "studyId",
    default: undefined,
});

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

            if(response?.data?.study) {
                // do some computations, e.g. generate arquero table of a received record list...
                const s: Study = {
                    ...response.data.study,
                    aStudyId: 2
                };
                return s;
            }
            return {
                aStudyId: 0
            } as Study;
        } else {
            return {
                aStudyId: 0
            } as Study;
        }
    }
});
