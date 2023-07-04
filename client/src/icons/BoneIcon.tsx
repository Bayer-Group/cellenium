import Icon from "./svg/bone.svg";

interface Icon {
  size: number;
}

const BoneIcon = ({ size }: Icon) => {
  return (
    <div style={{ width: size }}>
      <img src={Icon} alt="bone icon" />
    </div>
  );
};

export default BoneIcon;
