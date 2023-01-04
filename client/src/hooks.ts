import {useEffect, useMemo, useRef, useState} from 'react';
import * as aq from 'arquero';
import ColumnTable from 'arquero/dist/types/table/column-table';
import {useExpressionByOmicsIdsQuery} from "./generated/types";
import Table from "arquero/dist/types/table/table";

export function useExpressionValues() {

    const {data, loading} = useExpressionByOmicsIdsQuery({
        variables: {
            studyLayerId: 1,
            omicsIds: 116
        }
    });
    if (data?.expressionByOmicsIdsList) {
        let t = aq.from(data?.expressionByOmicsIdsList);
        // records contain two arrays, sampleIds and values, which have the same length and corresponding values.
        // unfold both arrays and put values side by side in new rows.
        const sampleIds = t.unroll('studySampleIds').select({studySampleIds: 'studySampleId'});
        t = t.unroll('values').select({values: 'value', omicsId: 'omicsId'});
        t = t.assign(sampleIds).reify()
        return {
            table: t,
            loading
        };
    }
    return {
        table: undefined,
        loading
    };
}