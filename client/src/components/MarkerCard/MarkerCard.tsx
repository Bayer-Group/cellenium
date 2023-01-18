import React from 'react';
import {DifferentialMarkerFragment} from "../../generated/types";
import {Badge, Card, Container, createStyles, Group, Spoiler, Stack, Text} from "@mantine/core";
import {InlineFoldChangePlot} from "../InlineFoldChangePlot/InlineFoldChangePlot";
import {useRecoilState} from "recoil";
import {
    annotationGroupIdState,
    pageState,
    selectedAnnotationState,
    selectedGenesState,
    userGenesState
} from "../../atoms";
import {useNavigate} from "react-router-dom";

interface Props {
    data: DifferentialMarkerFragment;
}

const useStyles = createStyles((theme) => ({
    main: {
        backgroundColor: "white",
        cursor: 'pointer',
        "&:hover": {
            backgroundColor: theme.colors.gray[1]
        }
    }
}));


const MarkerCard = ({data}: Props) => {
    const {classes, cx} = useStyles();
    const [page, setPage] = useRecoilState(pageState);
    const [selectedAnnotation, setSelectedAnnotation] = useRecoilState(selectedAnnotationState);
    const [annotationGroupId, setAnnotationGroupId] = useRecoilState(annotationGroupIdState);
    const [userGenes, setUserGenes] = useRecoilState(userGenesState)
    const [selectedGenes, setSelectedGenes] = useRecoilState(selectedGenesState);
    const navigate = useNavigate();

    function showGeneInStudy() {
        setPage('CellMarkerAnalysis');
        setUserGenes([data.omics]);
        setSelectedGenes([data.omics]);
        setAnnotationGroupId(data.annotationValue.annotationGroup.annotationGroupId)
        setSelectedAnnotation(data.annotationValueId);
        navigate(`/study/${data.study.studyId}`);
    }

    return (
        <Card pr={0} onClick={() => showGeneInStudy()} shadow="sm" p="lg" radius="md" withBorder>
            <Card.Section className={classes.main} withBorder inheritPadding py="xs">
                <Group position={'left'} spacing={'xs'} noWrap={true} pr={20}>
                    <Stack spacing={4} justify={'flex-end'} align={'flex-start'} style={{minWidth: 90}}>
                        <Text span fw={100} lineClamp={(1)} size={'xs'}>Cell
                            annotation</Text>
                        <Text fw={100} size={'xs'}>Study</Text>
                    </Stack>
                    <Stack spacing={4} justify={'flex-start'}>
                        <Text fw={700}
                              lineClamp={1}
                              size={'xs'}>{`${data.annotationValue.displayValue} (${data.annotationValue.annotationGroup.displayGroup})`}</Text>
                        <Text fw={700} lineClamp={1} size={'xs'}>{data.study.studyName}</Text>
                    </Stack>
                </Group>
            </Card.Section>
            <Stack pt={10} pb={10}>
                <Group position={'left'} spacing={'xs'} noWrap={true}>
                    <Stack align={'flex-start'}>
                        <Text fw={100} size={'xs'}>Gene</Text>
                        <Text size={'xs'} fw={100} lineClamp={1}>p-value (adj.)</Text>
                        <Text size={'xs'} fw={100}>log2FC</Text>
                    </Stack>
                    <Stack>
                        <Text size={'xs'} fw={700}>{data.omics.displaySymbol}</Text>
                        <Text size={'xs'} fw={700}>{data.pvalueAdj.toExponential(2)}</Text>
                        <Text size={'xs'} fw={700}>{data.log2Foldchange.toFixed(2)}</Text>
                    </Stack>
                    <Container>
                        <InlineFoldChangePlot studyId={data.study.studyId} annotationValueId={data.annotationValueId}
                                              pval={data.pvalueAdj}
                                              log2fc={data.log2Foldchange}/>
                    </Container>
                </Group>
            </Stack>

            <Card.Section withBorder inheritPadding py="xs">
                <Spoiler maxHeight={19} showLabel={"Show more"} hideLabel={'hide'} style={{
                    fontSize: 12,
                }}>
                    <Group position={'left'} spacing={3}>
                        {[1, 2, 3, 4].map((item) => {
                            let badges = [<Badge/>]
                            return badges;
                        })}
                    </Group>
                </Spoiler>

            </Card.Section>
        </Card>
    );
};

export {MarkerCard};
