import Icon from "./svg/breast.svg";
interface Icon {
  size: number;
}
const BreastIcon = ({ size }: Icon) => {
  return (
    <div style={{ width: size }}>
        <img src={Icon} alt="breast icon" />
    </div>
  );
};

export default BreastIcon;
