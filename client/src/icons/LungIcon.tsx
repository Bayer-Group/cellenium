import React from "react";
import { ReactComponent as Icon } from "./svg/lung.svg";
interface IMouseIcon {
  size: number;
}
const LungIcon = ({ size }: IMouseIcon) => {
  return (
    <div style={{ width: size }}>
      <Icon />
    </div>
  );
};

export default LungIcon;
