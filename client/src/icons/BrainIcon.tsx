import Icon from './svg/brain.svg';
interface Icon {
  size: number;
}
const BrainIcon = ({ size }: Icon) => {
  return (
    <div style={{ width: size }}>
      <img src={Icon} alt="brain icon" />
    </div>
  );
};

export default BrainIcon;
