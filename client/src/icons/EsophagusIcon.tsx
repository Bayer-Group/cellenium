import React from "react";
import { ReactComponent as Icon } from "./svg/esophagus.svg";
interface Icon {
  size: number;
}
const EsophagusIcon = ({ size }: Icon) => {
  return (
    <div style={{ width: size }}>
      <Icon />
    </div>
  );
};

export default EsophagusIcon;
