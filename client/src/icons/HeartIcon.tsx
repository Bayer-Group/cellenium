import React from "react";
import { ReactComponent as Heart } from "./svg/heart.svg";
interface IMouseIcon {
  size: number;
}
const HeartIcon = ({ size }: IMouseIcon) => {
  return (
    <div style={{ width: size }}>
      <Heart />
    </div>
  );
};

export default HeartIcon;
