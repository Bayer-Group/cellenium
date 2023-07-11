import Rat from "./svg/rat.svg";

interface Icon {
  size: number;
}

const RatIcon = ({ size }: Icon) => {
  return (
    <img src={Rat} className="h-auto" style={{ width: size }} alt="rat icon" />
  );
};

export default RatIcon;
