import {ActionIcon, TextInput, TextInputProps, useMantineTheme} from '@mantine/core';
import {IconPlus, IconSearch} from '@tabler/icons';

 function AddGene(props: TextInputProps) {
    const theme = useMantineTheme();

    return (
        <TextInput
            icon={<IconSearch size={18} stroke={1.5}/>}
            radius="md"
            size="xs"
            rightSection={
                <ActionIcon size={25} radius="md" color={theme.primaryColor} variant="filled">
                    <IconPlus size={14} stroke={1.5}/>
                </ActionIcon>
            }
            placeholder="Enter gene(s)"
            rightSectionWidth={30}
            {...props}
        />
    );
}
export {AddGene}