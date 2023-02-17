import React, {useEffect, useMemo, useState} from 'react';
import {Container, Grid, Loader, Space} from "@mantine/core";
import {StudyInfoFragment, useStudiesQuery} from "../../generated/types";
import {OntologyItem} from "../../model";
import {generateOntologyTrees} from "../../pages/helper";
import {SearchBar} from "./SearchBar";


function StudySearchBar({
                            onStudyListUpdate
                        }: {
    onStudyListUpdate: (studies: StudyInfoFragment[]) => void
}) {
    const {data, loading} = useStudiesQuery();
    const allStudies = useMemo(() => data?.studyOverviewsList && data.studyOverviewsList.map(study => ({
        ...study,
        allOntCodes: [
            study.studyOntologyList.map(ont => ont.ontCodes),
            study.studyOntologyList.map(ont => ont.parentIds)
        ].flat(2).filter(ontCode => !!ontCode)
    })), [data]);

    const ontologyTrees: Map<string, OntologyItem> | undefined = useMemo(
        () => data?.treeOntologiesList && generateOntologyTrees(data.treeOntologiesList),
        [data]);
    const [searchOntCodes, setSearchOntCodes] = useState<string[]>([]);

    const filteredStudies = useMemo(() => {
        if (!allStudies) {
            return undefined;
        }
        if (searchOntCodes.length === 0) {
            return allStudies;
        }
        return allStudies.filter(study => searchOntCodes.find(searchOntCode => study.allOntCodes.indexOf(searchOntCode) > -1));
    }, [data?.studyOverviewsList, searchOntCodes]);
    useEffect(() => onStudyListUpdate(filteredStudies || []), [filteredStudies]);

    if (loading) {
        return <Loader variant={'dots'}/>;
    }
    return <SearchBar ontologies={ontologyTrees} onSearchElementsUpdate={setSearchOntCodes}/>;
};

export default StudySearchBar;