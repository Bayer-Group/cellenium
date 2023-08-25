import Icon from './svg/lung.svg';
interface IMouseIcon {
  size: number;
}
const LungIcon = ({ size }: IMouseIcon) => {
  return (
    <div style={{ width: size }}>
      <img src={Icon} alt="lung icon" />
    </div>
  );
};

export default LungIcon;
