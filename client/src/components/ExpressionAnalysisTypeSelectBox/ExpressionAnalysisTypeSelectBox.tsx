import React from 'react';
import {Select} from "@mantine/core";
import {SelectBoxItem} from "../../model";
type Props = {
    handleSelection: any;
    selection: string;
    options: SelectBoxItem[];
}
const ExpressionAnalysisTypeSelectBox = ({selection, options, handleSelection}:Props) => {
    return (
        <Select value={selection} onChange={handleSelection} data={options} />
    );
};

export default ExpressionAnalysisTypeSelectBox;