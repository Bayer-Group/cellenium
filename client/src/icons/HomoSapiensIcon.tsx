import React from "react";
import { ReactComponent as HomoSapiens } from "./svg/homo_sapiens.svg";

interface IHomoSapiens {
  size: number;
}
const HomoSapiensIcon = ({ size }: IHomoSapiens) => {
  return (
      <HomoSapiens className="h-auto" style={{ width: size }}/>
  );
};

export default HomoSapiensIcon;
