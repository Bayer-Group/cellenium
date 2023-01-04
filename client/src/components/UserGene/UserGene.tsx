import React from 'react';
import {ActionIcon, CloseButton, Group, Text} from "@mantine/core";
import {IconDroplet} from "@tabler/icons";

const UserGene = ({gene}: any) => {
    return (
        <Group position={'apart'}>
            <Group>
                <CloseButton size={'xs'} iconSize={15}/>
                <Text size={'xs'}>{gene.display_symbol}</Text>
            </Group>
            {/* eslint-disable-next-line react/jsx-no-undef */}
            <ActionIcon variant="default" size={'xs'}>
                <IconDroplet/>
            </ActionIcon>
        </Group>
    );
};

export default UserGene;
