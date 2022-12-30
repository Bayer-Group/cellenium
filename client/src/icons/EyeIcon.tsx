import React from "react";
import { ReactComponent as Icon } from "./svg/eye.svg";
interface Icon {
  size: number;
}
const EyeIcon = ({ size }: Icon) => {
  return (
    <div style={{ width: size }}>
      <Icon />
    </div>
  );
};

export default EyeIcon;
