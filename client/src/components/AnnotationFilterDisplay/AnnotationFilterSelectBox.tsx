import { useMemo } from "react";
import { MultiSelect, SelectItem } from "@mantine/core";
import { useRecoilState, useRecoilValue } from "recoil";
import { selectedAnnotationFilterState, studyState } from "../../atoms";

const AnnotationFilterSelectBox = () => {
  const study = useRecoilValue(studyState);
  const [selectedAnnotationFilter, setSelectedAnnotationFilter] =
    useRecoilState(selectedAnnotationFilterState);
  const annotations: SelectItem[] = useMemo(() => {
    const anns: SelectItem[] = [];
    if (study) {
      study.annotationGroupMap.forEach((value) => {
        const groupName = value.displayGroup;
        value.annotationValuesList.map((av) =>
          anns.push({
            value: av.annotationValueId.toString(),
            label: av.displayValue,
            group: groupName,
          }),
        );
      });
    }
    return anns;
  }, [study]);

  return (
    <MultiSelect
      style={{ width: 210, maxWidth: 210, minWidth: 210 }}
      value={selectedAnnotationFilter.map((s) => s.toString())}
      onChange={(value) => {
        setSelectedAnnotationFilter(value.map((v) => parseInt(v)));
      }}
      searchable
      placeholder="Select"
      data={annotations}
    />
  );
};

export default AnnotationFilterSelectBox;
