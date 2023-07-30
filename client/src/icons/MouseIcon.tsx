import MusMusculus from './svg/mus_musculus.svg';

interface IMouseIcon {
  size: number;
}
const MouseIcon = ({ size }: IMouseIcon) => {
  return <img src={MusMusculus} className="h-auto" style={{ width: size }} alt="mus musculus icon" />;
};

export default MouseIcon;
