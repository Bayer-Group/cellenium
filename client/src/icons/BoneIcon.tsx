import React from "react";
import { ReactComponent as Icon } from "./svg/bone.svg";

interface Icon {
  size: number;
}

const BoneIcon = ({ size }: Icon) => {
  return (
    <div style={{ width: size }}>
      <Icon />
    </div>
  );
};

export default BoneIcon;
