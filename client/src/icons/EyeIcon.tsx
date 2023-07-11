import Icon from "./svg/eye.svg";
interface Icon {
  size: number;
}
const EyeIcon = ({ size }: Icon) => {
  return (
    <div style={{ width: size }}>
      <img src={Icon} alt="eye icon" />
    </div>
  );
};

export default EyeIcon;
