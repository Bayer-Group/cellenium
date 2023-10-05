import { ActionIcon, Anchor, Badge, Card, Grid, Group, Spoiler, Text } from '@mantine/core';
import { IconExternalLink } from '@tabler/icons-react';
import { Link } from 'react-router-dom';
import { StudyInfoFragment } from '../../generated/types';
import { ontology2Color } from '../../utils/helper';

function isValidEmail(email: string) {
  const emailRegex = /^[a-zA-Z0-9._%+-]+@[a-zA-Z0-9.-]+\.[a-zA-Z]{2,}$/;
  return emailRegex.test(email);
}

function isValidUrl(urlString: string) {
  try {
    return Boolean(new URL(urlString));
  } catch (e) {
    return false;
  }
}

function metadataKeyValues(study: StudyInfoFragment) {
  if (study.metadata) {
    return Object.keys(study.metadata)
      .map((key) => {
        const k = `${study.studyId}-${key}-label`;
        if (isValidUrl(study.metadata[key])) {
          return (
            <Badge key={k} color="gray">
              <Anchor color="gray" component={Link} to={study.metadata[key] as string}>
                {key}: {study.metadata[key]}
              </Anchor>
            </Badge>
          );
        }

        if (isValidEmail(study.metadata[key])) {
          return (
            <Badge key={k} color="gray">
              <Anchor color="gray" href={`mailto:${study.metadata[key]}`} target="_blank">
                {key}: {study.metadata[key]}
              </Anchor>
            </Badge>
          );
        }

        return (
          <Badge key={k} color="gray">
            {key}: {study.metadata[key]}
          </Badge>
        );
      })
      .sort();
  }
  return [];
}

export function StudyCard({ study, detailed }: { study: StudyInfoFragment; detailed?: boolean }) {
  const newStudyUrl = `study/${study.studyId}`;

  return (
    <Card shadow={detailed ? undefined : 'sm'} p="lg" radius="md" withBorder={!detailed}>
      <Card.Section withBorder inheritPadding py="xs">
        <Grid columns={12}>
          <Grid.Col span={8}>
            <Anchor component={Link} to={detailed ? '' : newStudyUrl} color="dark">
              <Text align="left" lineClamp={detailed ? undefined : 1} sx={detailed ? {} : { textOverflow: 'ellipsis', overflow: 'hidden' }} weight={800}>
                {study.studyName}
              </Text>
            </Anchor>
          </Grid.Col>
          <Grid.Col span={4}>
            <Group position="right">
              <Badge variant="light" color="gray">
                {Math.round(study.cellCount / 1000)}k cells
              </Badge>
              {/* eslint-disable-next-line react/jsx-no-undef */}
              {study.externalWebsite && (
                <ActionIcon
                  variant="subtle"
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
      <Text mt="sm" mb="sm" color="dimmed" size="sm" lineClamp={detailed ? undefined : 3} align="left">
        {study.description}
      </Text>
      <Card.Section withBorder inheritPadding py="xs">
        <Spoiler
          maxHeight={detailed ? 500 : 25}
          showLabel="Show more"
          hideLabel="hide"
          style={{
            fontSize: 12,
          }}
        >
          <Group position="left" spacing={3}>
            {study.studyOntologyList &&
              study.studyOntologyList.map((item) => {
                if (item.labels !== null) {
                  return item.labels.map((label) => (
                    <Badge color={ontology2Color(item.ontology)} key={`${study.studyId}-${label}-badge`}>
                      {label}
                    </Badge>
                  ));
                }
                return null;
              })}
            {metadataKeyValues(study)}
          </Group>
        </Spoiler>
      </Card.Section>
    </Card>
  );
}
