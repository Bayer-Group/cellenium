import React from "react";
import { ReactComponent as Icon } from "./svg/breast.svg";
interface Icon {
  size: number;
}
const BreastIcon = ({ size }: Icon) => {
  return (
    <div style={{ width: size }}>
      <Icon />
    </div>
  );
};

export default BreastIcon;
