import React from 'react';
import {ColorSwatch, createStyles, Grid, Text} from "@mantine/core";
import {useRecoilState} from "recoil";
import {highlightAnnotationState, selectedAnnotationState} from "../../atoms";


const useStyles = createStyles((theme) => ({
    hovered: {
        backgroundColor: theme.colors.gray[3],
        borderRadius: theme.radius.xs
    },
    clicked: {
        backgroundColor: theme.colors.blue[1],
        borderRadius: theme.radius.xs
    }
}));

type Props = {
    label: string;
    color: string;
    annotationId: number;
}

function Annotation({label, color, annotationId}: Props) {
    const {classes, cx} = useStyles();
    const [highlight, setHighlight] = useRecoilState(highlightAnnotationState);
    const [selected, setSelected] = useRecoilState(selectedAnnotationState);

    return (
        <Grid columns={12} pl={10} gutter={0} sx={{cursor: 'pointer'}} justify={'space-between'} align={'center'}
              onMouseOver={() => setHighlight(annotationId)}
              onClick={() => {
                  if (highlight === selected) {
                      setSelected(0)
                  } else {
                      setSelected(annotationId)
                  }
              }}
              className={cx({
                  [classes.hovered]: annotationId === highlight,
                  [classes.clicked]: annotationId === selected
              })}
        >
            <Grid.Col span={10}>
                <Text title={label} size={'md'} fw={selected === annotationId ? 800 : "md"} lineClamp={1}>
                    {label}
                </Text>
            </Grid.Col>
            <Grid.Col span={2} style={{justifyContent: 'center'}}>
                <ColorSwatch key={color} color={color} size={selected === annotationId ? 12 : 15}
                             style={{border: selected === annotationId ? '2px solid black' : ''}}/>
            </Grid.Col>
        </Grid>
    );


};

export {Annotation};