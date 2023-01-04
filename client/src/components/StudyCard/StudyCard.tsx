import React from 'react';
import {ActionIcon, Anchor, Badge, Card, Grid, Group, Spoiler, Text} from '@mantine/core';
import {IconExternalLink} from "@tabler/icons";
import { Study} from "../../generated/types";

type Props = {
    studyId: number;
    studyName: string;
    cellCount: number;
    description: string;
    tissueNcitIds: string[];
    diseaseMeshIds: string[];
}
const StudyCard = ({studyId, studyName, cellCount, description, tissueNcitIds, diseaseMeshIds}: Props) => {
    return (
        <Card shadow="sm" p="lg" radius="md" withBorder>
            <Card.Section withBorder inheritPadding py="xs">
                <Grid columns={12}>
                    <Grid.Col span={8}>
                        <Anchor href={'deg'} color={'dark'}>
                            <Text align='left' lineClamp={1} sx={{textOverflow: 'ellipsis', overflow: 'hidden'}}
                                  weight={800}>{studyName}</Text>
                        </Anchor>
                    </Grid.Col>
                    <Grid.Col span={4}>
                        <Group position={'right'}>
                            <Badge variant={'light'} color={'gray'}>{Math.round(cellCount / 1000)}k
                                cells</Badge>
                            {/* eslint-disable-next-line react/jsx-no-undef */}
                            <ActionIcon variant={'subtle'} onClick={() => {
                                window.open("www.google.de", "_blank")
                            }}>
                                <IconExternalLink/>
                            </ActionIcon>
                        </Group>
                    </Grid.Col>
                </Grid>
            </Card.Section>
            <Text mt="sm" mb='sm' color="dimmed" size="sm" lineClamp={3} align={'left'}>
                {description}
            </Text>
            <Card.Section withBorder inheritPadding py="xs">
                <Spoiler maxHeight={19} showLabel={"Show more"} hideLabel={'hide'} style={{
                    fontSize: 12,
                }}>
                    <Group position={'left'} spacing={3}>
                        {tissueNcitIds && tissueNcitIds.map((tissue: string) => {
                            return (<Badge size='sm' color="primary" variant="outline">
                                {tissue}
                            </Badge>)
                        })}
                        {diseaseMeshIds && diseaseMeshIds.map((disease: string) => {
                            if (disease.length > 0) return (<Badge size='sm' color="red" variant="outline">
                                {disease}
                            </Badge>)
                        })}
                    </Group>
                </Spoiler>

            </Card.Section>
        </Card>
    );
};

export {StudyCard};