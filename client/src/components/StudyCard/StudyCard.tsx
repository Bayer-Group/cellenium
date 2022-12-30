import React from 'react';
import {ActionIcon, Badge, Card, Grid, Group, Spoiler, Text} from '@mantine/core';
import {IconExternalLink} from "@tabler/icons";

const StudyCard = ({content}: any) => {
    return (
        <Card shadow="sm" p="lg" radius="md" withBorder>
            <Card.Section withBorder inheritPadding py="xs">
                <Grid columns={12}>
                    <Grid.Col span={8}>
                        <Text align='left' lineClamp={1} sx={{textOverflow: 'ellipsis', overflow: 'hidden'}}
                              weight={800}>{content.name}</Text>
                    </Grid.Col>
                    <Grid.Col span={4}>
                        <Group position={'right'}>
                            <Badge variant={'light'} color={'gray'}>{Math.round(content.cellCount / 1000)}k cells</Badge>
                            {/* eslint-disable-next-line react/jsx-no-undef */}
                            <ActionIcon variant={'subtle'} onClick={()=>{
                                window.open(content.source, "_blank")
                            }}>
                                <IconExternalLink/>
                            </ActionIcon>
                        </Group>
                    </Grid.Col>
                </Grid>
            </Card.Section>
            <Text mt="sm" mb='sm' color="dimmed" size="sm" lineClamp={3} align={'left'}>
                {content.description}
            </Text>
            <Card.Section withBorder inheritPadding py="xs">
                <Spoiler maxHeight={19} showLabel={"Show more"} hideLabel={'hide'} style={{
                    fontSize: 12,
                }}>
                    <Group position={'left'} spacing={3}>
                        {content.tissues.map((tissue: string) => {
                            return (<Badge size='sm' color="primary" variant="outline">
                                {tissue}
                            </Badge>)
                        })}
                        {content.diseases.map((disease: string) => {
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