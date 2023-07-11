import Icon from "./svg/ovary.svg";

interface Icon {
  size: number;
}

const OvaryIcon = ({ size }: Icon) => {
  return (
    <div style={{ width: size }}>
      <img src={Icon} alt="ovary icon" />
    </div>
  );
};

export default OvaryIcon;
