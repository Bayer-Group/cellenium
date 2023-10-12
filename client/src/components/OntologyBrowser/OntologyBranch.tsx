import { useCallback, useState } from 'react';
import { OntologyItem } from '../../model';
import { OntologyNode } from './OntologyNode';
import { StudyOverview } from '../../generated/types';

export function OntologyBranch({
  item,
  level,
  handleAddOntologyItem,
  studies,
}: {
  item: OntologyItem;
  level: number;
  handleAddOntologyItem: (item: OntologyItem) => void;
  studies?: StudyOverview & { allOntCodes: string[] }[];
}) {
  const [selected, setSelected] = useState(false);
  const hasChildren = !!(item && item.children && item.children.length > 0);

  const toggle = useCallback(() => {
    setSelected((prev) => !prev);
  }, [setSelected]);

  return (
    <>
      <OntologyNode
        studies={studies}
        item={item}
        selected={selected}
        hasChildren={hasChildren}
        level={level}
        onToggle={toggle}
        handleAddOntologyItem={handleAddOntologyItem}
      />
      {selected && hasChildren && item && item.children
        ? item.children.map((child) => {
            return <OntologyBranch studies={studies} handleAddOntologyItem={handleAddOntologyItem} key={child.unique_id} item={child} level={level + 1} />;
          })
        : null}
    </>
  );
}
