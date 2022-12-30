import React from "react";
import { ReactComponent as Icon } from "./svg/immune_system.svg";
interface Icon {
  size: number;
}
const ImmuneIcon = ({ size }: Icon) => {
  return (
    <div style={{ width: size }}>
      <Icon />
    </div>
  );
};

export default ImmuneIcon;
