import {useCallback, useState} from "react";
import {NavBar} from "../components";
import {
    Button,
    Checkbox,
    Container,
    Group,
    Loader,
    Modal,
    Select,
    Space,
    Stack,
    Text,
    Textarea,
    TextInput,
} from "@mantine/core";
import {useForm} from "@mantine/form";
import {IconDotsVertical, IconPlus, IconX} from "@tabler/icons-react";
import {
    InputMaybe,
    StudyAdminDetailsFragment,
    useCreateStudyUploadMutation,
    useStudyAdminListQuery,
    useStudyDeleteMutation,
    useStudyUpdateMutation,
} from "../generated/types";
import DataTable from "react-data-table-component";
import {showNotification} from "@mantine/notifications";
import {Form} from "react-router-dom";

function DeleteStudyModal({
                              opened,
                              reset,
                              study,
                          }: {
    opened: boolean;
    reset: () => void;
    study: StudyAdminDetailsFragment | undefined;
}) {
    const [confirmationInput, setConfirmationInput] = useState("");
    const [deleteStudyMutation, {loading}] = useStudyDeleteMutation();

    const deleteStudy = useCallback(() => {
        if (study === undefined) {
            return;
        }
        // TODO: delete s3 file
        deleteStudyMutation({
            variables: {
                studyId: study.studyId,
            },
        })
            .then(() => {
                setConfirmationInput("");
                reset();
            })
            .catch((reason) => {
                showNotification({
                    title: "Could not delete study",
                    message: reason.message,
                    color: "red",
                });
            });
    }, [study, deleteStudyMutation, reset]);

    return (
        <Modal opened={opened} onClose={reset}>
            <Stack>
                <Text weight="bold" size="xl">
                    Delete Study
                </Text>
                <Text>
                    Are you sure you want to delete the study <b>{study?.studyName}</b>?
                </Text>
                <Text>
                    To confirm deletion, type <i>permanently delete</i> in the text input
                    field.
                </Text>
                <TextInput
                    value={confirmationInput}
                    onChange={(event) => setConfirmationInput(event.currentTarget.value)}
                    placeholder="permanently delete"
                />
                <Group position="right">
                    <Button
                        color="red"
                        loading={loading}
                        disabled={confirmationInput !== "permanently delete" || loading}
                        onClick={deleteStudy}
                    >
                        Delete
                    </Button>
                </Group>
            </Stack>
        </Modal>
    );
}

function CreateStudyModal({
                              opened,
                              reset,
                          }: {
    opened: boolean;
    reset: () => void;
}) {
    const form = useForm({
        initialValues: {
            studyName: "",
            filetype: ".h5ad",
        },
        validate: (values) => {
            const errors: Record<string, string> = {};
            if (values.studyName === "") {
                errors.studyName = "Study name is required";
            }
            return errors;
        },
    });

    const [createStudyUploadMutation, {loading}] =
        useCreateStudyUploadMutation();

    const createStudy = useCallback(() => {
        createStudyUploadMutation({
            variables: {
                studyName: form.values.studyName,
                filetype: form.values.filetype,
            },
        })
            .then(() => {
                form.reset();
                reset();
            })
            .catch((reason) => {
                showNotification({
                    title: "Could not create study",
                    message: reason.message,
                    color: "red",
                });
            });
    }, [form, createStudyUploadMutation, reset]);

    return (
        <Modal opened={opened} onClose={reset}>
            <Stack>
                <Text weight="bold" size="xl">
                    Create Study
                </Text>
                <Form>
                    <TextInput label="Study Title" {...form.getInputProps("studyName")} />
                    <Select
                        data={[".h5ad", ".h5mu"]}
                        label="Filetype"
                        placeholder="Select a filetype"
                        {...form.getInputProps("filetype")}
                    />
                </Form>
                <Group position="right">
                    <Button
                        color="blue"
                        onClick={createStudy}
                        loading={loading}
                        disabled={!form.isValid()}
                    >
                        Create
                    </Button>
                </Group>
            </Stack>
        </Modal>
    );
}

export function NewStudyAdmin() {
    const [newStudyModalOpen, setNewStudyModalOpen] = useState(false);
    // const createTempCredentials = () => {
    //   createS3TempCredentialsMutation()
    //     .then((r) => {
    //       if (r?.data?.createS3TempCredentials?.strings) {
    //         setS3TempCredentials(r.data.createS3TempCredentials.strings);
    //         setNewStudyModalOpen(true);
    //       }
    //     })
    //     .catch((reason) => {
    //       showNotification({
    //         title: "Could not create temporary credentials",
    //         message: reason.message,
    //         color: "red",
    //       });
    //     });
    // };

    const {data, loading, refetch} = useStudyAdminListQuery();
    const [selectedStudy, setSelectedStudy] = useState<
        StudyAdminDetailsFragment | undefined
    >(undefined);
    const [selectedDeleteStudy, setSelectedDeleteStudy] = useState<
        StudyAdminDetailsFragment | undefined
    >(undefined);

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
            externalWebsite: selectedRow?.externalWebsite || "",
        });
    };

    const [studyUpdateMutation, {loading: studyUpdateLoading}] =
        useStudyUpdateMutation();
    const submit = () => {
        if (selectedStudy) {
            const splitToArray = (s: string) => {
                const a = s.split(";").map((s) => s.trim());
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
                    readerPermissions: splitToArray(
                        form.values.readerPermissions,
                    ) as InputMaybe<string[]>,
                    adminPermissions: splitToArray(
                        form.values.adminPermissions,
                    ) as InputMaybe<string[]>,
                    tissueNcitIds: splitToArray(form.values.tissueNcitIds) as InputMaybe<
                        string[]
                    >,
                    diseaseMeshIds: splitToArray(
                        form.values.diseaseMeshIds,
                    ) as InputMaybe<string[]>,
                    visible: form.values.visible,
                    externalWebsite: form.values.externalWebsite,
                },
            })
                .then(() => {
                    selectStudy(undefined);
                    void refetch();
                })
                .catch((reason) => {
                    showNotification({
                        title: "Could not save study changes",
                        message: reason.message,
                        color: "red",
                    });
                });
        }
    };

    const resetDeleteModal = useCallback(() => {
        setSelectedDeleteStudy(undefined);
        void refetch();
    }, []);

    const form = useForm({
        initialValues: {
            studyName: "",
            description: "",
            readerPermissions: "",
            adminPermissions: "",
            tissueNcitIds: "",
            diseaseMeshIds: "",
            visible: false,
            externalWebsite: "",
        },
        validate: {},
    });

    const columns = [
        {
            name: "ID",
            selector: (row: StudyAdminDetailsFragment) => row.studyId,
            sortable: true,
        },
        {
            name: "Title",
            selector: (row: StudyAdminDetailsFragment) => row.studyName,
            sortable: true,
        },
        {
            name: "Filename",
            selector: (row: StudyAdminDetailsFragment) => row.filename,
            sortable: true,
        },
        {
            name: "Your Role",
            selector: (row: StudyAdminDetailsFragment) =>
                row.adminPermissionGranted
                    ? "Admin"
                    : row.readerPermissionGranted
                        ? "View"
                        : "No Access",
            sortable: true,
        },
        {
            name: "Import Status",
            selector: (row: StudyAdminDetailsFragment) =>
                row.importStarted
                    ? row.importFailed
                        ? "Failed"
                        : row.importFinished
                            ? "Imported"
                            : "Unknown Error"
                    : row.importFinished
                        ? "Imported"
                        : row.importFailed ? "Failed":"Not Started",
            sortable: true,
        },
        {
            name: "",
            cell: (row: StudyAdminDetailsFragment) =>
                !row.adminPermissionGranted ? null : (
                    <IconDotsVertical
                        onClick={() => selectStudy(row)}
                        style={{cursor: "pointer"}}
                    />
                ),
        },
        {
            name: "",
            cell: (row: StudyAdminDetailsFragment) =>
                !row.adminPermissionGranted ? null : (
                    <IconX
                        onClick={() => setSelectedDeleteStudy(row)}
                        style={{cursor: "pointer"}}
                    />
                ),
        },
    ];

    return (
        <Container fluid={true}>
            <Modal
                opened={selectedStudy !== undefined}
                onClose={() => {
                    setSelectedStudy(undefined);
                }}
                size="80vw"
            >
                <Stack>
                    <Text weight="bold" size="xl">
                        Edit Study
                    </Text>
                    <Form>
                        <TextInput label="Title" {...form.getInputProps("studyName")} />
                        <Textarea
                            label="Description"
                            {...form.getInputProps("description")}
                        />
                        <TextInput
                            label="Reader Permissions, separate multiple groups / usernames with ;"
                            {...form.getInputProps("readerPermissions")}
                        />
                        <TextInput
                            label="Admin Permissions, separate multiple groups / usernames with ;"
                            {...form.getInputProps("adminPermissions")}
                        />
                        <Checkbox
                            mt="md"
                            label="Study is visible"
                            {...form.getInputProps("visible", {type: "checkbox"})}
                        />
                        <TextInput
                            label="Tissue NCIT IDs, separate multiple with ;"
                            {...form.getInputProps("tissueNcitIds")}
                        />
                        <TextInput
                            label="Disease MeSH IDs, separate multiple with ; and use the pseudo-ID HEALTHY to indicate 'not diseased'"
                            {...form.getInputProps("diseaseMeshIds")}
                        />
                        <TextInput
                            label="External Website"
                            {...form.getInputProps("externalWebsite")}
                        />
                        <Group position="right" mt="md">
                            <Button
                                disabled={!selectedStudy?.adminPermissionGranted}
                                onClick={submit}
                                loading={studyUpdateLoading}
                            >
                                Save Changes
                            </Button>
                        </Group>
                    </Form>
                </Stack>
            </Modal>
            <CreateStudyModal
                opened={newStudyModalOpen}
                reset={() => {
                    setNewStudyModalOpen(false);
                    void refetch();
                }}
            />
            <NavBar/>
            <Space h="xl"/>
            <Stack px="md">
                <Group position="right">
                    {data?.userStudyUploadConfigured && (
                        <Button onClick={() => setNewStudyModalOpen(true)}>
                            <Group spacing="xs">
                                <IconPlus/>
                                <span>New Study</span>
                            </Group>
                        </Button>
                    )}
                </Group>

                <DeleteStudyModal
                    study={selectedDeleteStudy}
                    reset={resetDeleteModal}
                    opened={selectedDeleteStudy !== undefined}
                />

                {loading && <Loader variant={"dots"} color={"gray"} size={"xl"}/>}
                {!loading ? <DataTable
                    data={data?.studyAdminDetailsList || []}
                    columns={columns}
                    onSelectedRowsChange={(state) =>
                        selectStudy(
                            state.selectedRows.length === 1
                                ? state.selectedRows[0]
                                : undefined,
                        )
                    }
                /> : null}
            </Stack>
        </Container>
    );
}
