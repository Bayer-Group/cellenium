import React from 'react';
import {DifferentialMarkerFragment} from "../../generated/types";
import {Badge, Card, Container, createStyles, Group, Spoiler, Stack, Text} from "@mantine/core";
import {InlineFoldChangePlot} from "../InlineFoldChangePlot/InlineFoldChangePlot";

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

    return (
        <Card pr={0} onClick={() => alert(1)} shadow="sm" p="lg" radius="md" withBorder>
            <Card.Section className={classes.main} withBorder inheritPadding py="xs">
                <Group position={'apart'} spacing={'xs'} noWrap={true} pr={20}>
                    <Stack spacing={4}>
                        <Badge variant="outline">cell
                            annotation</Badge>
                        <Badge variant="outline">study</Badge>
                    </Stack>
                    <Stack spacing={4}>
                        <Text fw={700} lineClamp={1}>{data.annotationValue.displayValue}</Text>
                        <Text fw={700} lineClamp={1}>{data.study.studyName}</Text>
                    </Stack>
                </Group>
            </Card.Section>
            <Stack pt={10} pb={10}>
                <Group position={'left'} spacing={'xs'} noWrap={true}>
                    <Stack>
                        <Text fw={700} size={'sm'}>Gene</Text>
                        <Text size={'sm'} fw={700}>p-value (adjusted)</Text>
                        <Text size={'sm'} fw={700}>log2FC</Text>
                    </Stack>
                    <Stack>
                        <Text size={'sm'}>{data.omics.displaySymbol}</Text>
                        <Text size={'sm'}>{data.pvalueAdj.toExponential(2)}</Text>
                        <Text size={'sm'}>{data.log2Foldchange.toFixed(2)}</Text>
                    </Stack>
                    <Container>
                        <InlineFoldChangePlot/>
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
