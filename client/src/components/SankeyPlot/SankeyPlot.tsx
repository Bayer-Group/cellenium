import { useMemo } from 'react';
import { DefaultLink, DefaultNode, ResponsiveSankey } from '@nivo/sankey';
import { Center, Loader, MantineColor, Stack } from '@mantine/core';
import { StudyAnnotationFrontendValue, useAnnotationValueCoocurrenceQuery } from '../../generated/types';

export function SankeyPlot({
  annotationGroupId1,
  annotationGroupId2,
  studyId,
  annotationValues1,
  annotationValues2,
}: {
  annotationGroupId1: number;
  annotationGroupId2: number;
  annotationValues1: StudyAnnotationFrontendValue[];
  annotationValues2: StudyAnnotationFrontendValue[];
  studyId: number;
}) {
  const { data, loading } = useAnnotationValueCoocurrenceQuery({
    variables: {
      studyId,
      annotationGroupId1,
      annotationGroupId2,
    },
  });
  const sankeyData = useMemo(() => {
    if (data) {
      const nodes: {
        id: number;
        nodeColor: MantineColor;
        label: string;
      }[] = [];
      const vals1: number[] = [];
      const vals2: number[] = [];
      annotationValues1.forEach((nd) => {
        nodes.push({
          id: nd.annotationValueId,
          nodeColor: nd.color,
          label: nd.displayValue,
        });
        vals1.push(nd.annotationValueId);
      });
      annotationValues2.forEach((nd) => {
        nodes.push({
          id: nd.annotationValueId,
          nodeColor: nd.color,
          label: nd.displayValue,
        });
        vals2.push(nd.annotationValueId);
      });
      const links: {
        source: number;
        target: number;
        value: number;
      }[] = [];
      data.annotationValueCoocurrenceList
        .filter((nd) => {
          return vals1.includes(nd.annotationValueId1) && vals2.includes(nd.annotationValueId2);
        })
        .forEach((nd) => {
          links.push({
            source: nd.annotationValueId1,
            target: nd.annotationValueId2,
            value: nd.occurrence,
          });
        });
      return {
        nodes,
        links,
      };
    }
    return {
      nodes: [],
      links: [],
    };
  }, [data, annotationValues1, annotationValues2]);

  return (
    <Center w="100%" h="100%">
      {loading && <Loader variant="dots" color="blue" size="xl" />}
      {sankeyData && sankeyData.nodes.length > 0 && (
        <Stack w="90%" h="90%">
          <ResponsiveSankey
            data={sankeyData as unknown as { nodes: DefaultNode[]; links: DefaultLink[] }}
            label="label"
            margin={{ top: 40, right: 250, bottom: 40, left: 250 }}
            align="justify"
            colors={{ scheme: 'category10' }}
            nodeOpacity={1}
            nodeHoverOthersOpacity={0.35}
            nodeThickness={18}
            nodeSpacing={24}
            nodeBorderWidth={0}
            nodeBorderColor={{
              from: 'color',
              modifiers: [['darker', 0.8]],
            }}
            nodeBorderRadius={3}
            linkOpacity={0.5}
            linkHoverOthersOpacity={0.1}
            linkContract={3}
            enableLinkGradient
            labelPosition="outside"
            labelOrientation="horizontal"
            labelPadding={16}
            labelTextColor={{
              from: 'color',
              modifiers: [['darker', 1]],
            }}
          />
        </Stack>
      )}
    </Center>
  );
}
