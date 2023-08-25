import { ActionIcon, Anchor, Badge, Card, Grid, Group, Spoiler, Text } from '@mantine/core';
import { IconExternalLink } from '@tabler/icons-react';
import { Link } from 'react-router-dom';
import { StudyInfoFragment } from '../../generated/types';
import { ontology2Color } from '../../pages/helper';

const StudyCard = ({ study }: { study: StudyInfoFragment }) => {
  const newStudyUrl = `study/${study.studyId}`;

  return (
    <Card shadow="sm" p="lg" radius="md" withBorder>
      <Card.Section withBorder inheritPadding py="xs">
        <Grid columns={12}>
          <Grid.Col span={8}>
            <Anchor component={Link} to={newStudyUrl} color={'dark'}>
              <Text align="left" lineClamp={1} sx={{ textOverflow: 'ellipsis', overflow: 'hidden' }} weight={800}>
                {study.studyName}
              </Text>
            </Anchor>
          </Grid.Col>
          <Grid.Col span={4}>
            <Group position={'right'}>
              <Badge variant={'light'} color={'gray'}>
                {Math.round(study.cellCount / 1000)}k cells
              </Badge>
              {/* eslint-disable-next-line react/jsx-no-undef */}
              {study.externalWebsite && (
                <ActionIcon
                  variant={'subtle'}
                  onClick={() => {
                    window.open(study.externalWebsite, '_blank');
                  }}
                >
                  <IconExternalLink />
                </ActionIcon>
              )}
            </Group>
          </Grid.Col>
        </Grid>
      </Card.Section>
      <Text mt="sm" mb="sm" color="dimmed" size="sm" lineClamp={3} align={'left'}>
        {study.description}
      </Text>
      <Card.Section withBorder inheritPadding py="xs">
        <Spoiler
          maxHeight={25}
          showLabel={'Show more'}
          hideLabel={'hide'}
          style={{
            fontSize: 12,
          }}
        >
          <Group position={'left'} spacing={3}>
            {study.studyOntologyList &&
              study.studyOntologyList.map((item) => {
                if (item.labels !== null) {
                  let badges = item.labels.map((label) => <Badge color={ontology2Color(item.ontology)}>{label}</Badge>);
                  return badges;
                }
                return null;
              })}
          </Group>
        </Spoiler>
      </Card.Section>
    </Card>
  );
};

export { StudyCard };
