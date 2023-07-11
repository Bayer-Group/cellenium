import { ActionIcon, Loader, Stack, Text } from "@mantine/core";
import { IconPlus } from "@tabler/icons-react";
import memoize from "memoize-one";
import { useRecoilState } from "recoil";
import {
  selectedGenesState,
  userGenesState,
  userGeneStoreCounterColor,
  userGeneStoreOpenState,
} from "../../atoms";
import { Omics } from "../../model";
import _ from "lodash";
import { useCorrelatedgenesQuery } from "../../generated/types";
import DataTable from "react-data-table-component";

const customStyles = {
  table: {
    style: {
      backgroundColor: "transparent",
      marginRight: "10",
      overflow: "hidden",
    },
  },
  header: {
    style: {
      paddingLeft: "2px",
      backgroundColor: "transparent",
    },
  },
  head: {
    style: {
      paddingLeft: "2px",
      backgroundColor: "transparent",
    },
  },
  headRow: {
    style: {
      paddingLeft: "2px",
      backgroundColor: "transparent",
    },
  },
  rows: {
    style: {
      minHeight: "72px", // override the row height
      backgroundColor: "transparent",
    },
  },
  headCells: {
    style: {
      paddingLeft: "2px", // override the cell padding for head cells
      paddingRight: "2px",
      backgroundColor: "transparent",
    },
  },
  cells: {
    style: {
      paddingLeft: "2px", // override the cell padding for data cells
      paddingRight: "2px",
      backgroundColor: "transparent",
    },
  },
};
const columns = memoize((clickHandler) => [
  {
    name: "gene",
    selector: (row: any) => row.displaySymbol,
    sortable: true,
    width: "80px",
  },
  {
    name: "r",
    selector: (row: any) => row.r.toFixed(2),
    sortable: true,
    width: "70px",
  },
  {
    name: "",
    cell: (row: any) => {
      let gene = {
        omicsId: row.omicsId,
        displayName: row.displayName,
        displaySymbol: row.displaySymbol,
        omicsType: row.omicsType,
        value: row.displaySymbol,
      };
      return (
        <ActionIcon
          color={"blue.3"}
          onClick={() => clickHandler(gene)}
          size="xs"
          variant={"default"}
        >
          <IconPlus size={12} color={"black"} />
        </ActionIcon>
      );
    },
    width: "20px",
  },
]);

type Props = {
  omicsId: number;
  studyId: number;
};

const CorrelationTable = ({ omicsId, studyId }: Props) => {
  const [userGenes, setUserGenes] = useRecoilState(userGenesState);
  const [selectedGenes, setSelectedGenes] = useRecoilState(selectedGenesState);
  const [, setIndicatorColor] = useRecoilState(userGeneStoreCounterColor);
  // const study = useRecoilValue(studyState);
  const [, setStoreOpen] = useRecoilState(userGeneStoreOpenState);
  const { data, loading } = useCorrelatedgenesQuery({
    variables: {
      omicsId: omicsId,
      studyId: studyId,
    },
  });

  function handleClick(gene: Omics) {
    let check = userGenes.filter((g) => g.omicsId === gene.omicsId);
    if (check.length === 0) {
      setSelectedGenes([...selectedGenes, gene]);

      setIndicatorColor("pink");
      setUserGenes(_.union(userGenes, [gene]));
      setStoreOpen(false);
      setTimeout(() => {
        setIndicatorColor("blue");
      }, 200);
    }
  }

  if (loading) {
    return (
      <Stack align={"center"}>
        <Text color={"gray"} size={"xs"}>
          A genome-wide correlation analysis takes some seconds. Please stay
          patient! We are working on a speed-up in the meantime!
        </Text>
        <Loader variant={"dots"} color={"gray"} />
      </Stack>
    );
  }

  return (
    <Stack justify={"flex-start"} align={"flex-start"}>
      {data && data.getCorrelatedGenesList.length > 0 && (
        <DataTable
          dense
          columns={columns(handleClick)}
          data={data.getCorrelatedGenesList}
          defaultSortFieldId={3}
          defaultSortAsc={false}
          customStyles={customStyles}
          fixedHeader
          fixedHeaderScrollHeight="100%"
          noDataComponent={<Text>No data.</Text>}
        />
      )}
      {data && data.getCorrelatedGenesList.length === 0 && (
        <Text color={"dimmed"} size={"xs"}>
          No correlated genes with Pearson r&gt;=0.2 found.
        </Text>
      )}
    </Stack>
  );
};

export { CorrelationTable };
