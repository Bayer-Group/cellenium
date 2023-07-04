import {ActionIcon, Autocomplete, AutocompleteItem} from '@mantine/core';
import {useEffect, useState} from "react";
import {IconX} from "@tabler/icons-react";
import {useRecoilValue} from "recoil";
import {studyState,} from "../../atoms";
import * as aq from 'arquero';
import {Omics} from "../../model";

interface Props {
    selection: Omics | null;
    onSelectionChange: (omics: Omics | null) => void;
}

function SingleGeneSelection({selection, onSelectionChange}: Props) {
    const [offerings, setOfferings] = useState<Omics[]>([]);
    const [value, setValue] = useState(selection ? selection.displaySymbol : '');
    useEffect(() => {
        if (value !== selection?.displaySymbol) {
            if (selection && value === '') {
                setValue(selection.displaySymbol);
            }
        }
    }, [selection]);
    const study = useRecoilValue(studyState);

    function handleChange(inputString: string) {
        let newOfferings: Omics[] = [];
        if (inputString.length > 0) {
            // @ts-ignore
            newOfferings = study?.studyOmicsTable.filter(aq.escape(t => aq.op.includes(t.displaySymbol.toLowerCase(), inputString.toLowerCase(), 0))).objects();
        }
        setOfferings(newOfferings)
        setValue(inputString)
    }

    function handleItemSubmit(item: Omics) {
        setOfferings([])
        if (item) {
            onSelectionChange(item);
        }
    }

    return (
        <Autocomplete
            value={value}
            radius={0}
            size={'xs'}
            onChange={handleChange}
            onBlur={() => handleItemSubmit(offerings[0])}
            data={offerings as AutocompleteItem[]}
            onItemSubmit={(item: any) => {
                handleItemSubmit(item)
            }}
            rightSection={
                <ActionIcon onClick={() => {
                    setValue('')
                    setOfferings([])
                    onSelectionChange(null);
                }
                }>
                    <IconX/>
                </ActionIcon>
            }
        />
    );
}

export {SingleGeneSelection}
