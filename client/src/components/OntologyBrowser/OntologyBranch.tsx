import { useState } from 'react';
import { OntologyItem } from '../../model';
import OntologyNode from './OntologyNode';

const OntologyBranch = ({ item, level, handleAddOntologyItem }: { item: OntologyItem; level: number; handleAddOntologyItem: Function }) => {
  const [selected, setSelected] = useState(false);
  const hasChildren = item && item.children && item.children.length > 0 ? true : false;

  const renderChildren = () => {
    if (hasChildren) {
      return (
        item &&
        item.children &&
        item.children.map((child) => {
          return <OntologyBranch handleAddOntologyItem={handleAddOntologyItem} key={child.unique_id} item={child} level={level + 1} />;
        })
      );
    }
    return null;
  };
  const toggle = () => {
    setSelected((prev) => !prev);
  };
  return (
    <>
      <OntologyNode item={item} selected={selected} hasChildren={hasChildren} level={level} onToggle={toggle} handleAddOntologyItem={handleAddOntologyItem} />
      {selected && renderChildren()}
    </>
  );
};

export default OntologyBranch;
