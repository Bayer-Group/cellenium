import {useEffect, useMemo} from 'react';
import * as aq from 'arquero';
import {InputMaybe, ProjectionType, useExpressionByOmicsIdsQuery} from "./generated/types";
import {ExpressionTable} from "./model";
import {useRecoilState, useRecoilValue, useResetRecoilState, useSetRecoilState} from "recoil";
import {
    annotationGroupIdState,
    highlightAnnotationState,
    pageState,
    selectedAnnotationState,
    selectedGenesState,
    studyIdState,
    studyLayerIdState,
    studyState,
    userGenesState
} from "./atoms";
import {useLocation, useNavigate, useParams} from "react-router-dom";

export function useExpressionValues(omicsIds: number[], subsampling: boolean) {
    const studyLayerId = useRecoilValue(studyLayerIdState);

    const {data, loading} = useExpressionByOmicsIdsQuery({
        variables: {
            studyLayerId,
            omicsIds,
            subsamplingProjection: (subsampling ? ProjectionType.Umap : null ) as InputMaybe<ProjectionType>
        },
        skip: omicsIds.length === 0
    });
    return useMemo(() => {
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
    }, [data, loading]);
}

export function useSetStudyFromUrl() {
    const {
        studyId: studyIdUrlParam
    } = useParams<{
        studyId: string
    }>();
    const queryParams = new URLSearchParams(useLocation().search);
    const navigate = useNavigate();
    const page = queryParams.get('page');

    const setValidParam = (queryParam: string, setter: (x: any) => void, defaultValue?: any) => {
        const value = queryParams.get(queryParam);
        if (value) {
            if (queryParam.endsWith('Id')) {
                const numericValue = Number(value);
                setter(numericValue);
            } else {
                setter(value);
            }
        } else if (defaultValue) {
            setter(defaultValue);
        }
    };

    const studyIdUrlParamInt = studyIdUrlParam && parseInt(studyIdUrlParam);
    const [studyId, setStudyId] = useRecoilState(studyIdState);
    const resetStatesOnStudyChange = [
        useResetRecoilState(selectedGenesState),
        useResetRecoilState(selectedAnnotationState),
        useResetRecoilState(highlightAnnotationState)
    ];

    useEffect(() => {
        if (studyIdUrlParamInt) {
            if (studyIdUrlParamInt > 0 && studyIdUrlParamInt !== studyId) {
                setStudyId(studyIdUrlParamInt);
                resetStatesOnStudyChange.forEach(atomReset => atomReset());
            }
        }
    }, [studyIdUrlParam, studyId]);

    // after the study is loaded, we set valid defaults or URL query params:
    const study = useRecoilValue(studyState);
    const setPage = useSetRecoilState(pageState);
    const setAnnotationGroupId = useSetRecoilState(annotationGroupIdState);
    const setUserGenes = useSetRecoilState(userGenesState);
    const setSelectedGenes = useSetRecoilState(selectedGenesState);
    const setSelectedAnnotation = useSetRecoilState(selectedAnnotationState);
    useEffect(() => {
        if (study && study.studyId === studyIdUrlParamInt) {
            setValidParam('page', setPage);
            setValidParam('annotationGroupId', setAnnotationGroupId, study.annotationGroupsList[0].annotationGroupId);
            setValidParam('annotationValueId', setSelectedAnnotation);
            setValidParam('omicsId', (omicsId: number) => {
                const o = study.studyOmicsMap.get(omicsId);
                setUserGenes(o ? [o] : []);
                setSelectedGenes(o ? [o] : []);
            });
            // clear query parameters, as we don't plan to update them and it's not nice to leave stale data in the URL
            navigate(`/study/${studyId}`, {replace: true});
        }
    }, [study]);
}

