import React from "react";
import { ReactComponent as Icon } from "./svg/BinaryTree.svg";

interface Icon {
  size?: string;
  className?: string
}

const BinaryTreeIcon = ({ size, className }: Icon) => {
  return (
      <Icon className={className}/>
  );
};

export default BinaryTreeIcon;
