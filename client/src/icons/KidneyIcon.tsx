import Kidneys from "./svg/kidneys.svg";
interface Icon {
  size: number;
}
const KidneyIcon = ({ size }: Icon) => {
  return (
    <div style={{ width: size }}>
      <img src={Kidneys} alt="kidney icon" />
    </div>
  );
};

export default KidneyIcon;
