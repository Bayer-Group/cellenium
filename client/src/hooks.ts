import { useEffect, useMemo } from 'react';
import * as aq from 'arquero';
import { SetterOrUpdater, useRecoilState, useRecoilValue, useResetRecoilState, useSetRecoilState } from 'recoil';
import { useLocation, useNavigate, useParams } from 'react-router-dom';
import { InputMaybe, useExpressionByOmicsIdsQuery } from './generated/types';
import { ExpressionTable } from './model';
import {
  annotationGroupIdState,
  annotationSecondaryGroupIdState,
  highlightAnnotationState,
  pageState,
  selectedAnnotationState,
  selectedGenesState,
  selectedProjectionState,
  studyIdState,
  studyLayerIdState,
  studyState,
  userGenesState,
} from './atoms';

export function useExpressionValues(omicsIds: number[], subsampling: boolean) {
  const studyLayerId = useRecoilValue(studyLayerIdState);
  const projection = useRecoilValue(selectedProjectionState);

  const { data, loading } = useExpressionByOmicsIdsQuery({
    variables: {
      studyLayerId,
      omicsIds,
      subsamplingProjection: (subsampling ? projection.split(':').pop() : null) as InputMaybe<string>,
    },
    skip: omicsIds.length === 0,
  });
  return useMemo(() => {
    if (data?.expressionByOmicsIdsList) {
      let t = aq.from(data?.expressionByOmicsIdsList);
      // records contain two arrays, sampleIds and values, which have the same length and corresponding values.
      // unfold both arrays and put values side by side in new rows.
      if (t.numRows() > 0) {
        const sampleIds = t.unroll('studySampleIds').select({ studySampleIds: 'studySampleId' });
        t = t.unroll('values').select({ values: 'value', omicsId: 'omicsId' });
        t = t.assign(sampleIds).reify();
        return {
          table: ExpressionTable.definedTable(t),
          loading,
        };
      }
    }
    return {
      table: undefined,
      loading,
    };
  }, [data, loading]);
}

export function useSetStudyFromUrl() {
  const { studyId: studyIdUrlParam } = useParams<{
    studyId: string;
  }>();
  const queryParams = new URLSearchParams(useLocation().search);
  const navigate = useNavigate();
  // const page = queryParams.get('page');

  const setValidParam = <T>(queryParam: string, setter: SetterOrUpdater<T> | ((v: T) => void), defaultValue?: T) => {
    const value = queryParams.get(queryParam);
    if (value) {
      if (queryParam.endsWith('Id')) {
        const numericValue = Number(value);
        setter(numericValue as T);
      } else {
        setter(value as T);
      }
    } else if (defaultValue) {
      setter(defaultValue);
    }
  };

  const studyIdUrlParamInt = studyIdUrlParam && parseInt(studyIdUrlParam, 10);
  const [studyId, setStudyId] = useRecoilState(studyIdState);
  const resetStatesOnStudyChange = [
    useResetRecoilState(selectedGenesState),
    useResetRecoilState(selectedAnnotationState),
    useResetRecoilState(highlightAnnotationState),
  ];

  useEffect(() => {
    if (studyIdUrlParamInt) {
      if (studyIdUrlParamInt > 0 && studyIdUrlParamInt !== studyId) {
        setStudyId(studyIdUrlParamInt);
        resetStatesOnStudyChange.forEach((atomReset) => atomReset());
      }
    }
  }, [studyIdUrlParam, studyId]);

  // after the study is loaded, we set valid defaults or URL query params:
  const study = useRecoilValue(studyState);
  const setPage = useSetRecoilState(pageState);
  const setAnnotationGroupId = useSetRecoilState(annotationGroupIdState);
  const setSecondaryAnnotationGroupId = useSetRecoilState(annotationSecondaryGroupIdState);
  const setUserGenes = useSetRecoilState(userGenesState);
  const setSelectedGenes = useSetRecoilState(selectedGenesState);
  const setSelectedAnnotation = useSetRecoilState(selectedAnnotationState);
  const setSelectedProjection = useSetRecoilState(selectedProjectionState);
  useEffect(() => {
    if (study && study.studyId === studyIdUrlParamInt) {
      setValidParam('page', setPage);
      setValidParam('annotationGroupId', setAnnotationGroupId, study.annotationGroupsList[0].annotationGroupId);
      setValidParam('annotationValueId', setSelectedAnnotation);
      setSecondaryAnnotationGroupId(undefined);
      setValidParam('omicsId', (omicsId: number) => {
        const o = study.studyOmicsMap.get(omicsId);
        setUserGenes(o ? [o] : []);
        setSelectedGenes(o ? [o] : []);
      });
      setValidParam('projection', setSelectedProjection, study.projections[0]);
      // clear query parameters, as we don't plan to update them and it's not nice to leave stale data in the URL
      navigate(`/study/${studyId}`, { replace: true });
    }
  }, [study]);
}
