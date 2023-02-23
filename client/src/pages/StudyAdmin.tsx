import React, {useEffect, useMemo, useState} from 'react';
import {NavBar} from "../components";
import {
    ActionIcon,
    Checkbox,
    Container,
    Grid,
    Group,
    Loader,
    Space,
    TextInput,
    Button,
    Box,
    Divider
} from "@mantine/core";
import {useForm} from '@mantine/form';
import {
    InputMaybe,
    StudyAdminDetailsFragment,
    useStudyAdminListQuery, useStudyUpdateMutation
} from "../generated/types";
import DataTable from 'react-data-table-component';
import {showNotification} from "@mantine/notifications";

function StudyAdmin() {
    const {data, loading, refetch} = useStudyAdminListQuery();

    const columns = [
        {
            name: 'ID',
            selector: (row: StudyAdminDetailsFragment) => row.studyId,
            sortable: true,
        },
        {
            name: 'Title',
            selector: (row: StudyAdminDetailsFragment) => row.studyName,
            sortable: true,
        },
        {
            name: 'Filename',
            selector: (row: StudyAdminDetailsFragment) => row.filename,
            sortable: true,
        },
        {
            name: 'Your Role',
            selector: (row: StudyAdminDetailsFragment) => row.adminPermissionGranted ? "Admin" : (row.readerPermissionGranted ? "View" : "No Access"),
            sortable: true,
        },
    ];

    const [selectedStudy, setSelectedStudy] = useState<StudyAdminDetailsFragment | undefined>(undefined);

    const form = useForm({
        initialValues: {
            studyName: '',
            description: '',
            readerPermissions: '',
            adminPermissions: '',
            tissueNcitIds: '',
            diseaseMeshIds: '',
            visible: false,
            externalWebsite: '',
        },
        validate: {},
    });


    const selectStudy = (selectedRow: StudyAdminDetailsFragment | undefined) => {
        setSelectedStudy(selectedRow);
        form.setValues({
            studyName: selectedRow?.studyName || "",
            description: selectedRow?.description || "",
            readerPermissions: (selectedRow?.readerPermissions || []).join("; "),
            adminPermissions: (selectedRow?.adminPermissions || []).join("; "),
            tissueNcitIds: (selectedRow?.tissueNcitIds || []).join("; "),
            diseaseMeshIds: (selectedRow?.diseaseMeshIds || []).join("; "),
            visible: selectedRow?.visible || false,
            externalWebsite: selectedRow?.externalWebsite || ""
        });
    }


    const [studyUpdateMutation, {
        loading: studyUpdateLoading
    }] = useStudyUpdateMutation();
    const submit = () => {
        console.log('SUBMIT')
        if (selectedStudy) {
            const splitToArray = (s: string) => {
                const a = s.split(";").map(s => s.trim());
                if (a.length === 1 && a[0] === "") {
                    return null;
                }
                return a;
            }

            studyUpdateMutation({
                variables: {
                    studyId: selectedStudy.studyId,
                    studyName: form.values.studyName,
                    description: form.values.description,
                    readerPermissions: splitToArray(form.values.readerPermissions) as InputMaybe<string[]>,
                    adminPermissions: splitToArray(form.values.adminPermissions) as InputMaybe<string[]>,
                    tissueNcitIds: splitToArray(form.values.tissueNcitIds) as InputMaybe<string[]>,
                    diseaseMeshIds: splitToArray(form.values.diseaseMeshIds) as InputMaybe<string[]>,
                    visible: form.values.visible,
                    externalWebsite: form.values.externalWebsite
                }
            }).then(() => {
                selectStudy(undefined);
                refetch();
            }).catch(reason => {
                showNotification({
                    title: 'Could not save study changes',
                    message: reason.message,
                    color: 'red'
                });
            });
        }
    };

    return <Container fluid={true}>
        <NavBar/>
        <Space h="xl"/>

        {loading && <Loader variant={'dots'} color={'gray'} size={'xl'}/>}
        <DataTable data={data?.studyAdminDetailsList || []}
                   columns={columns}
                   selectableRows
                   selectableRowsSingle
                   onSelectedRowsChange={state => selectStudy(state.selectedRows.length === 1 ? state.selectedRows[0] : undefined)}
        />

        <Space h="xl"/>
        <Box>
            <form>
                <TextInput label="Title"
                           {...form.getInputProps('studyName')}
                />
                <TextInput label="Description"
                           {...form.getInputProps('description')}
                />
                <TextInput label="Reader Permissions, separate multiple groups / usernames with ;"
                           {...form.getInputProps('readerPermissions')}
                />
                <TextInput label="Admin Permissions, separate multiple groups / usernames with ;"
                           {...form.getInputProps('adminPermissions')}
                />
                <Checkbox
                    mt="md"
                    label="Study is visible"
                    {...form.getInputProps('visible', {type: 'checkbox'})}
                />
                <TextInput label="Tissue NCIT IDs, separate multiple with ;"
                           {...form.getInputProps('tissueNcitIds')}
                />
                <TextInput
                    label="Disease MeSH IDs, separate multiple with ; and use the pseudo-ID HEALTHY to indicate 'not diseased'"
                    {...form.getInputProps('diseaseMeshIds')}
                />
                <TextInput label="External Website"
                           {...form.getInputProps('externalWebsite')}
                />
                <Group position="right" mt="md">
                    <Button disabled={!(selectedStudy?.adminPermissionGranted)} onClick={submit}
                            loading={studyUpdateLoading}>Save Changes</Button>
                </Group>
            </form>
        </Box>
    </Container>;
}

export default StudyAdmin;
