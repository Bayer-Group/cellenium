import React from "react";
import { ReactComponent as Fly } from "./svg/fly.svg";

interface Icon {
  size: number;
}
const FlyIcon = ({ size }: Icon) => {
  return (
      <Fly className="h-auto" style={{ width: size }} />
  );
};

export default FlyIcon;
