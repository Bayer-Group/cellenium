import React from "react";
import { ReactComponent as Icon } from "./svg/ovary.svg";
interface Icon {
  size: number;
}
const OvaryIcon = ({ size }: Icon) => {
  return (
    <div style={{ width: size }}>
      <Icon />
    </div>
  );
};

export default OvaryIcon;
