import React, { forwardRef, useMemo, useState } from "react";
import { Group, Select, Text } from "@mantine/core";
import { useRecoilState, useRecoilValue } from "recoil";
import {
  annotationGroupIdState,
  selectedAnnotationState,
  studyState,
} from "../../atoms";
import { IconCalculator } from "@tabler/icons-react";

interface ItemProps extends React.ComponentPropsWithoutRef<"div"> {
  value: string;
  label: string;
  differentialExpressionCalculated: boolean;
}

const SelectItem = forwardRef<HTMLDivElement, ItemProps>(
  ({ label, differentialExpressionCalculated, ...others }: ItemProps, ref) => (
    <div ref={ref} {...others}>
      <Group position={"apart"} align={"center"} noWrap>
        <Text>{label}</Text>
        {differentialExpressionCalculated && (
          <Group title="DEG calculated" align={"center"}>
            <IconCalculator color={"gray"} size={20} />
          </Group>
        )}
      </Group>
    </div>
  ),
);

function AnnotationGroupSelectBox() {
  const study = useRecoilValue(studyState);
  const [annotationGroupId, setAnnotationGroupId] = useRecoilState(
    annotationGroupIdState,
  );
  const [, setSelectedAnnotation] = useRecoilState(selectedAnnotationState);
  const [value, setValue] = useState<string | undefined>();
  const annotations: ItemProps[] = useMemo(() => {
    const anns: ItemProps[] = [];
    if (study) {
      study.annotationGroupMap.forEach((value, key) => {
        anns.push({
          value: key.toString(),
          label: value.displayGroup,
          differentialExpressionCalculated:
            value.differentialExpressionCalculated,
        });
      });
    }
    setValue(annotationGroupId?.toString());
    return anns;
  }, [study, annotationGroupId]);

  function update(value: string | null) {
    if (value) {
      setValue(value);
      setAnnotationGroupId(parseInt(value));
      setSelectedAnnotation(0);
    }
  }

  return (
    <Select
      style={{ maxWidth: 210, width: 210, minWidth: 210 }}
      value={value}
      onChange={(value) => update(value)}
      label="Select annotation group"
      labelProps={{ size: "xs" }}
      itemComponent={SelectItem}
      placeholder="Pick one"
      transitionProps={{
        duration: 80,
        timingFunction: "ease",
      }}
      data={annotations as any}
    />
  );
}

export { AnnotationGroupSelectBox };
