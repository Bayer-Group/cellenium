import React from "react";
import { ReactComponent as Rat } from "./svg/rat.svg";

interface Icon {
  size: number;
}
const RatIcon = ({ size }: Icon) => {
  return (
      <Rat className="h-auto" style={{ width: size }}/>
  );
};

export default RatIcon;
