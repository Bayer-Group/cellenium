import {ActionIcon, TextInput, TextInputProps, useMantineTheme} from '@mantine/core';
import {IconArrowLeft, IconArrowRight, IconSearch} from '@tabler/icons';

function SearchBar(props: TextInputProps) {
    const theme = useMantineTheme();

    return (
        <TextInput
            variant='default'
            label="filter studies by disease, tissue, species, title/description"
            styles={{
                label: {fontWeight:100, fontSize: '0.8rem', display: 'inline-block'}
            }}
            icon={<IconSearch size={18} stroke={1.5}/>}

            size="md"
            rightSection={
                <ActionIcon size={32} radius="sm" color={theme.primaryColor} variant="filled">
                    {theme.dir === 'ltr' ? (
                        <IconArrowRight size={18} stroke={1.5}/>
                    ) : (
                        <IconArrowLeft size={18} stroke={1.5}/>
                    )}
                </ActionIcon>
            }
            placeholder='lung, "multiple myelome", heart, mouse'
            rightSectionWidth={42}
            {...props}
        />
    );
}

export {SearchBar};