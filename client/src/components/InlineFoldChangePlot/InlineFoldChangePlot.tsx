import Plot from 'react-plotly.js';
import { Loader } from '@mantine/core';
import { useHalfAVolcanoQuery } from '../../generated/types';

export function InlineFoldChangePlot({
  studyId,
  annotationValueId,
  pval,
  log2fc,
}: {
  studyId: number;
  annotationValueId: number;
  pval: number;
  log2fc: number;
}) {
  const { data, loading } = useHalfAVolcanoQuery({
    variables: {
      studyId,
      annotationValueId,
    },
  });
  if (loading) return <Loader variant="dots" color="blue" />;
  // p values of 0 get set to 10e-100 to avoid infinity in the log10 scaled plot
  const pvalues = data?.differentialExpressionsList && data.differentialExpressionsList.map((e) => (e.pvalueAdj === 0 ? 10 ** -100 : e.pvalueAdj));
  const minPvalueAdj = (pvalues && Math.min(...pvalues)) || 10 ** -100;

  return (
    <Plot
      config={{
        modeBarButtons: false,
        displaylogo: false,
        staticPlot: true,
      }}
      data={[
        {
          x: data?.differentialExpressionsList.map((e) => e.log2Foldchange),
          y: pvalues?.map((e) => {
            const ppval = -1 * Math.log10(e);
            return ppval || 0;
          }),
          type: 'scatter',
          mode: 'markers',
          marker: { color: 'lightblue', size: 3 },
          showlegend: false,
        },
        {
          x: [log2fc],
          y: [-1 * Math.log10(pval || minPvalueAdj)],
          marker: {
            color: 'red',
            size: 5,
          },
          showlegend: false,
        },
      ]}
      layout={{
        width: 160,
        height: 100,
        margin: { l: 26, r: 5, t: 5, b: 20 },

        xaxis: {
          visible: true,
          fixedrange: true,
          showgrid: false,
          title: 'log2FC',
          titlefont: {
            size: 10,
            color: 'grey',
          },
          tickfont: {
            size: 6,
          },
        },
        yaxis: {
          visible: true,
          fixedrange: true,
          showgrid: false,
          title: '-log(padj)',
          titlefont: {
            family: 'Arial, sans-serif',
            size: 10,
            color: 'grey',
          },
          tickfont: {
            size: 6,
          },
        },
      }}
    />
  );
}
