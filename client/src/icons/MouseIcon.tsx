import React from "react";
import { ReactComponent as MusMusculus } from "./svg/mus_musculus.svg";

interface IMouseIcon {
  size: number;
}
const MouseIcon = ({ size }: IMouseIcon) => {
  return (
      <MusMusculus className="h-auto" style={{ width: size }} />
  );
};

export default MouseIcon;
