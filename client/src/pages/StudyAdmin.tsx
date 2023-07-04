import {useState} from 'react';
import {NavBar} from "../components";
import {Box, Button, Checkbox, Container, Group, Loader, Modal, Space, Stack, TextInput} from "@mantine/core";
import {useForm} from '@mantine/form';
import {
    InputMaybe,
    StudyAdminDetailsFragment,
    useCreateS3TempCredentialsMutation,
    useStudyAdminListQuery,
    useStudyUpdateMutation
} from "../generated/types";
import DataTable from 'react-data-table-component';
import {showNotification} from "@mantine/notifications";

function StudyAdmin() {
    const [tempCredentialsModalOpened, setTempCredentialsModalOpened] = useState(false);
    const [s3TempCredentials, setS3TempCredentials] = useState<string[]>([]);
    const [createS3TempCredentialsMutation, {
        loading: createS3TempCredentialsLoading
    }] = useCreateS3TempCredentialsMutation();
    const createTempCredentials = () => {
        createS3TempCredentialsMutation().then((r) => {
            if (r?.data?.createS3TempCredentials?.strings) {
                setS3TempCredentials(r.data.createS3TempCredentials.strings);
                setTempCredentialsModalOpened(true);
            }
        }).catch(reason => {
            showNotification({
                title: 'Could not create temporary credentials',
                message: reason.message,
                color: 'red'
            });
        });
    };


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
        if (selectedStudy) {
            const splitToArray = (s: string) => {
                const a = s.split(";").map(s => s.trim());
                if (a.length === 1 && a[0] === "") {
                    return null;
                }
                return a;
            };

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
        <Modal
            opened={tempCredentialsModalOpened}
            onClose={() => setTempCredentialsModalOpened(false)}
            title="AWS credentials for study upload"
            size={"50em"}
        >
            <Stack>
                <p>
                    An h5ad scRNA expression data file needs to fulfill certain requirements to be of use in cellenium.
                    Please use "h5ad_preparation.py" to ensure your file meets the requirements.
                </p>
                <p>
                    The credentials below are personalized to {s3TempCredentials[3]} and valid for 12 hours.
                    Access to the S3 bucket is logged. Please upload an h5ad study file as indicated below. After
                    successful file transfer, wait and refresh the study list to see the results.
                </p>
                <pre style={{overflow: "auto", width: "100%"}}>
{`export AWS_ACCESS_KEY_ID="${s3TempCredentials[0]}"
export AWS_SECRET_ACCESS_KEY="${s3TempCredentials[1]}"
export AWS_SESSION_TOKEN="${s3TempCredentials[2]}"

aws s3 ls ${s3TempCredentials[4]}
aws s3 cp --acl bucket-owner-full-control   new_study.h5ad   ${s3TempCredentials[4]}`}
           </pre>
            </Stack>
        </Modal>


        <NavBar/>
        <Space h="xl"/>

        {data?.userStudyUploadConfigured && <Group>
            <Button onClick={createTempCredentials} loading={createS3TempCredentialsLoading}>
                Study Upload: Create S3 Credentials
            </Button>
        </Group>}

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
