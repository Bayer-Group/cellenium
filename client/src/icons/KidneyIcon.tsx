import React from "react";
import { ReactComponent as Kidneys } from "./svg/kidneys.svg";
interface Icon {
  size: number;
}
const KidneyIcon = ({ size }: Icon) => {
  return (
    <div style={{ width: size }}>
      <Kidneys />
    </div>
  );
};

export default KidneyIcon;
