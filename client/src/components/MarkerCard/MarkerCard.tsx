import React from 'react';
import {DifferentialMarkerFragment} from "../../generated/types";
import {Card, Title} from "@mantine/core";

interface Props {
    data: DifferentialMarkerFragment;
}

const MarkerCard = ({data}: Props) => {
    return (
        <Card>
            <Title>{data.annotationValue.displayValue}</Title>
        </Card>
    );
};

export {MarkerCard};