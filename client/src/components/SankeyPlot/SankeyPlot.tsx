import { useMemo } from 'react';
import { useAnnotationValueCoocurrenceQuery } from '../../generated/types';
import { ResponsiveSankey } from '@nivo/sankey';
import { Center, Loader } from '@mantine/core';

interface Props {
  annotationGroupId1: number;
  annotationGroupId2: number;
  annotationValues1: any;
  annotationValues2: any;
  studyId: number;
}

const SankeyPlot = ({ annotationGroupId1, annotationGroupId2, studyId, annotationValues1, annotationValues2 }: Props) => {
  const { data, loading } = useAnnotationValueCoocurrenceQuery({
    variables: {
      studyId,
      annotationGroupId1,
      annotationGroupId2,
    },
  });
  const sankeyData = useMemo(() => {
    if (data) {
      const nodes: any = [];
      const vals1: number[] = [];
      const vals2: number[] = [];
      annotationValues1.map((nd: any) => {
        nodes.push({
          id: nd.annotationValueId,
          nodeColor: nd.color,
          label: nd.displayValue,
        });
        vals1.push(nd.annotationValueId);
      });
      annotationValues2.map((nd: any) => {
        nodes.push({
          id: nd.annotationValueId,
          nodeColor: nd.color,
          label: nd.displayValue,
        });
        vals2.push(nd.annotationValueId);
      });
      const links: any = [];
      data.annotationValueCoocurrenceList
        .filter((nd) => {
          return vals1.includes(nd.annotationValueId1) && vals2.includes(nd.annotationValueId2);
        })
        .map((nd) => {
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
  }, [data, annotationGroupId1, annotationGroupId2, annotationValues1, annotationValues2]);

  return (
    <Center style={{ width: '100%', height: '100%' }}>
      {loading && <Loader variant={'dots'} color={'gray'} size={'xl'} />}

      {sankeyData && sankeyData.nodes.length > 0 && (
        <div style={{ width: '90%', height: '100%' }}>
          <ResponsiveSankey
            data={sankeyData}
            label={'label'}
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
            enableLinkGradient={true}
            labelPosition="outside"
            labelOrientation="horizontal"
            labelPadding={16}
            labelTextColor={{
              from: 'color',
              modifiers: [['darker', 1]],
            }}
          />
        </div>
      )}
    </Center>
  );
};

export { SankeyPlot };

/*
 */
