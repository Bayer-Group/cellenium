import {useEffect, useMemo, useRef, useState} from 'react';
import * as aq from 'arquero';
import {useExpressionByOmicsIdsQuery} from "./generated/types";
import {ExpressionTable} from "./model";
import {useRecoilState, useRecoilValue, useResetRecoilState, useSetRecoilState} from "recoil";
import {
    highlightAnnotationState,
    selectedAnnotationState,
    selectedGenesState,
    studyIdState,
    studyLayerIdState
} from "./atoms";
import {useParams} from "react-router-dom";

export function useExpressionValues(omicsIds: number[]) {
    const studyLayerId = useRecoilValue(studyLayerIdState);

    const {data, loading} = useExpressionByOmicsIdsQuery({
        variables: {
            studyLayerId,
            omicsIds
        },
        skip: omicsIds.length === 0
    });
    if (data?.expressionByOmicsIdsList) {
        let t = aq.from(data?.expressionByOmicsIdsList);
        // records contain two arrays, sampleIds and values, which have the same length and corresponding values.
        // unfold both arrays and put values side by side in new rows.
        if (t.numRows() > 0) {
            const sampleIds = t.unroll('studySampleIds').select({studySampleIds: 'studySampleId'});
            t = t.unroll('values').select({values: 'value', omicsId: 'omicsId'});
            t = t.assign(sampleIds).reify()
            return {
                table: ExpressionTable.definedTable(t),
                loading
            };
        }
    }
    return {
        table: undefined,
        loading
    };
}

export function useSetStudyFromUrl() {
    const {studyId: studyIdUrlParam} = useParams<{ studyId: string }>();
    const [studyId, setStudyId] = useRecoilState(studyIdState);
    const resetStatesOnStudyChange = [
        useResetRecoilState(selectedGenesState),
        useResetRecoilState(selectedAnnotationState),
        useResetRecoilState(highlightAnnotationState)
    ];

    useEffect(() => {
        if (studyIdUrlParam) {
            const studyIdUrlParamInt = parseInt(studyIdUrlParam);
            if (studyIdUrlParamInt > 0 && studyIdUrlParamInt !== studyId) {
                setStudyId(studyIdUrlParamInt);
                resetStatesOnStudyChange.forEach(atomReset => atomReset());
            }
        }
    }, [studyIdUrlParam, studyId]);
}