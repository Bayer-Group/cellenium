import { Select, Text } from "@mantine/core";
import { useRecoilState, useRecoilValue } from "recoil";
import { selectedProjectionState, studyState } from "../../atoms";

const ProjectionSelectBox = () => {
  const study = useRecoilValue(studyState);
  const [projection, setProjection] = useRecoilState(selectedProjectionState);
  const niceLabel = (value: string) =>
    value
      .replace("umap", "UMAP")
      .replace("tsne", "t-SNE")
      .replace("pca", "PCA");
  const options = Array.from(study?.projections || []).map((p) => ({
    value: p,
    label: niceLabel(p),
  }));

  if (options.length === 1) {
    return <Text size={"xs"}>Projection: {options[0].label}</Text>;
  }

  return (
    <Select
      style={{ minWidth: 210, maxWidth: 210, width: 210 }}
      labelProps={{ size: "xs" }}
      label={"Select projection"}
      value={projection}
      onChange={(p) => setProjection(p || "")}
      data={options}
    />
  );
};

export default ProjectionSelectBox;
