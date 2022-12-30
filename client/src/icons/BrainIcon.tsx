import React from "react";
import { ReactComponent as Icon } from "./svg/brain.svg";
interface Icon {
  size: number;
}
const BrainIcon = ({ size }: Icon) => {
  return (
    <div style={{ width: size }}>
      <Icon />
    </div>
  );
};

export default BrainIcon;
