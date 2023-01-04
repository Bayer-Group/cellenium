import {useState} from 'react';
import {createStyles, Navbar, Title, Tooltip, UnstyledButton} from '@mantine/core';
import {ReactComponent as StudyAnalysisIcon} from "../../icons/study_analysis.svg";
import {ReactComponent as ExpressionAnalysisIcon} from "../../icons/expression_analysis.svg";
import {ReactComponent as CoExpressionAnalysisIcon} from "../../icons/coexpression_analysis.svg";
import {ReactComponent as CompareAnnotationsIcon} from "../../icons/annotation_comparison.svg";
import {ReactComponent as UserAnnotationIcon} from "../../icons/user_annotation.svg";

import {ReactComponent as CelleniumLogo} from "../../images/logo.svg";

const useStyles = createStyles((theme) => ({
    wrapper: {
        display: 'flex',
    },

    main: {
        flex: 1,
        borderLeft: '1px solid #e9efef',
        padding: theme.spacing.md,
        backgroundColor: theme.colorScheme === 'dark' ? theme.colors.dark[6] : theme.colors.gray[0],
    },
}));

type Props = {
    children?: JSX.Element[]|JSX.Element;
}

function RightSidePanel({children}:Props) {
    const {classes, cx} = useStyles();

    return (
        <Navbar height={'100vh'} width={{sm: 300}}>
            <Navbar.Section grow className={classes.wrapper}>
                <div className={classes.main}>
                    {children}
                </div>
            </Navbar.Section>
        </Navbar>
    );
}

export {RightSidePanel}